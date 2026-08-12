! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_solver_from_host
  ! Tests the radiative transfer solver that a host application updates
  !
  ! Without an argument, the program runs the standard tests. With an
  ! assertion code as the argument, it runs the failure test for that code.
  ! ``solver_from_host.sh`` drives the failure tests.

  use musica_assert,                   only : assert, die
  use musica_constants,                only : dk => musica_dk
  use musica_mpi,                      only : musica_mpi_finalize,            &
                                              musica_mpi_init

  implicit none

  ! shape of the radiation field that the test configurations describe
  integer, parameter :: n_interfaces = 5
  integer, parameter :: n_bins = 6
  ! contents of the test configurations
  integer, parameter :: n_reactions   = 2
  integer, parameter :: n_dose_rates  = 2
  ! the test configurations set the extraterrestrial flux to this value in
  ! every wavelength bin
  real(dk), parameter :: etfl = 1.0e14_dk
  ! cross section times quantum yield for each photolysis reaction
  real(dk), parameter :: xsqy( n_reactions ) = (/ 2.0_dk * 0.5_dk,            &
                                                  4.0_dk * 0.25_dk /)
  ! mid-point of each wavelength bin [nm]
  real(dk), parameter :: lambda( n_bins ) = (/ 425.0_dk, 475.0_dk, 525.0_dk,  &
                                               575.0_dk, 625.0_dk, 675.0_dk /)
  ! the two notch filters in the test configuration. The first passes every
  ! bin. The second passes the bins whose mid-point is above 500 nm.
  real(dk), parameter :: weight( n_bins, n_dose_rates ) =                     &
      reshape( (/ 1.0_dk, 1.0_dk, 1.0_dk, 1.0_dk, 1.0_dk, 1.0_dk,             &
                  0.0_dk, 0.0_dk, 1.0_dk, 1.0_dk, 1.0_dk, 1.0_dk /),          &
               (/ n_bins, n_dose_rates /) )
  real(dk), parameter :: earth_sun_distance = 0.9_dk
  real(dk), parameter :: tol = 1.0e-10_dk

  character(len=256) :: failure_test_code

  call musica_mpi_init( )
  if( command_argument_count( ) == 1 ) then
    ! A failure test stops inside the assertion that it triggers, so it never
    ! reaches the call to die below
    call get_command_argument( 1, failure_test_code )
    call failure_test( failure_test_code )
    call die( 197432558 )
  end if
  call assert( 806415392, command_argument_count( ) == 0 )
  call test_solver_from_host_t( )
  call test_core_with_host_radiation_field( )
  call test_updater_for_another_solver( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_solver_from_host_t( )
    ! Tests the solver in isolation, including the MPI transfer

    use musica_assert,                 only : assert, die
    use musica_config,                 only : config_t
    use musica_mpi,                    only : musica_mpi_bcast,               &
                                              musica_mpi_rank, MPI_COMM_WORLD
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_solver,                   only : solver_t, radiation_field_t
    use tuvx_solver_factory,           only : solver_allocate, solver_builder,&
                                              solver_type_name
    use tuvx_solver_from_host,         only : radiation_field_updater_t,      &
                                              solver_from_host_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t
    use tuvx_test_utils,               only : check_values

    character(len=*), parameter :: Iam = "from host solver tests"
    integer, parameter :: comm = MPI_COMM_WORLD

    character, allocatable :: buffer(:)
    class(grid_warehouse_t),     pointer :: grids
    class(profile_warehouse_t),  pointer :: profiles
    class(radiator_warehouse_t), pointer :: radiators
    class(solver_t),             pointer :: solver
    type(config_t) :: config, sub_config
    type(radiation_field_t), pointer :: field
    type(radiation_field_updater_t)  :: updater
    type(spherical_geometry_t), pointer :: geometry
    type(string_t) :: type_name
    integer :: i_interface, i_bin, pos, pack_size
    real(dk) :: fdr( n_interfaces, n_bins ), fup( n_interfaces, n_bins )
    real(dk) :: fdn( n_interfaces, n_bins ), edr( n_interfaces, n_bins )
    real(dk) :: eup( n_interfaces, n_bins ), edn( n_interfaces, n_bins )
    real(dk) :: zeros( n_interfaces, n_bins )

    call config%from_file( "test/data/solver_from_host.json" )
    call config%get( "grids", sub_config, Iam )
    grids => grid_warehouse_t( sub_config )
    call config%get( "profiles", sub_config, Iam )
    profiles => profile_warehouse_t( sub_config, grids )
    radiators => radiator_warehouse_t( )
    geometry => spherical_geometry_t( grids )

    ! build the solver on the primary rank and transfer it to the others
    if( musica_mpi_rank( comm ) == 0 ) then
      call sub_config%empty( )
      call sub_config%add( "type", "from host", Iam )
      solver => solver_builder( sub_config, grids, profiles )
      type_name = solver_type_name( solver )
      call assert( 316530884, type_name == "solver_from_host_t" )
      pack_size = type_name%pack_size( comm ) + solver%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack( buffer, pos, comm )
      call solver%mpi_pack( buffer, pos, comm )
      call assert( 493457630, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      solver => solver_allocate( type_name )
      call solver%mpi_unpack( buffer, pos, comm )
      call assert( 105833723, pos <= pack_size )
    end if
    deallocate( buffer )

    select type( solver )
    class is( solver_from_host_t )
      call assert( 282760469,                                                 &
                   solver%number_of_vertical_interfaces( ) == n_interfaces )
      call assert( 730128315, solver%number_of_wavelength_bins( ) == n_bins )
      updater = radiation_field_updater_t( solver )
    class default
      call die( 559971411 )
    end select

    zeros(:,:) = 0.0_dk

    ! a solver that the host application has never updated holds a radiation
    ! field of zero
    field => solver%update_radiation_field( 30.0_dk, n_interfaces - 1,        &
                                            geometry, grids, profiles,        &
                                            radiators )
    call check_values( 349517206, field%fdr_, zeros, tol )
    call check_values( 861885551, field%eup_, zeros, tol )
    deallocate( field )

    do i_bin = 1, n_bins
      do i_interface = 1, n_interfaces
        fdr( i_interface, i_bin ) = 0.1_dk * i_interface + 0.01_dk * i_bin
        fup( i_interface, i_bin ) = 0.001_dk * i_interface
        fdn( i_interface, i_bin ) = 0.002_dk * i_bin
        edr( i_interface, i_bin ) = 1.1_dk * i_interface
        eup( i_interface, i_bin ) = 2.2_dk * i_bin
        edn( i_interface, i_bin ) = 3.3_dk * i_interface * i_bin
      end do
    end do

    ! the actinic flux components alone; the irradiances default to zero
    call updater%update( direct_actinic_flux   = fdr,                         &
                         upward_actinic_flux   = fup,                         &
                         downward_actinic_flux = fdn )
    field => solver%update_radiation_field( 30.0_dk, n_interfaces - 1,        &
                                            geometry, grids, profiles,        &
                                            radiators )
    call check_values( 977055161, field%fdr_, fdr,   tol )
    call check_values( 306898257, field%fup_, fup,   tol )
    call check_values( 754266103, field%fdn_, fdn,   tol )
    call check_values( 584109199, field%edr_, zeros, tol )
    call check_values( 196485446, field%eup_, zeros, tol )
    call check_values( 926336041, field%edn_, zeros, tol )
    deallocate( field )

    ! all six components
    call updater%update( direct_actinic_flux   = fdr,                         &
                         upward_actinic_flux   = fup,                         &
                         downward_actinic_flux = fdn,                         &
                         direct_irradiance     = edr,                         &
                         upward_irradiance     = eup,                         &
                         downward_irradiance   = edn )
    field => solver%update_radiation_field( 30.0_dk, n_interfaces - 1,        &
                                            geometry, grids, profiles,        &
                                            radiators )
    call check_values( 133261945, field%fdr_, fdr, tol )
    call check_values( 645630290, field%fup_, fup, tol )
    call check_values( 475473386, field%fdn_, fdn, tol )
    call check_values( 987841731, field%edr_, edr, tol )
    call check_values( 817684827, field%eup_, eup, tol )
    call check_values( 365052674, field%edn_, edn, tol )

    ! the solver returns a copy, so the caller can free it safely
    deallocate( field )
    field => solver%update_radiation_field( 30.0_dk, n_interfaces - 1,        &
                                            geometry, grids, profiles,        &
                                            radiators )
    call check_values( 812420520, field%fdr_, fdr, tol )
    deallocate( field )

    deallocate( solver    )
    deallocate( geometry  )
    deallocate( radiators )
    deallocate( profiles  )
    deallocate( grids     )

  end subroutine test_solver_from_host_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_core_with_host_radiation_field( )
    ! Tests that TUV-x calculates photolysis rate constants and dose rates
    ! from a radiation field that the host application supplies

    use musica_assert,                 only : assert
    use musica_string,                 only : string_t
    use tuvx_constants,                only : hc
    use tuvx_core,                     only : core_t
    use tuvx_solver,                   only : radiation_field_t
    use tuvx_solver_from_host,         only : radiation_field_updater_t
    use tuvx_test_utils,               only : check_values

    class(core_t), pointer :: core
    type(radiation_field_updater_t) :: updater
    type(radiation_field_t) :: field
    type(string_t) :: config_path
    type(string_t), allocatable :: labels(:)
    logical :: found
    integer :: i_interface, i_bin, i_rate
    real(dk) :: fdr( n_interfaces, n_bins ), fup( n_interfaces, n_bins )
    real(dk) :: fdn( n_interfaces, n_bins )
    real(dk) :: edr( n_interfaces, n_bins ), eup( n_interfaces, n_bins )
    real(dk) :: edn( n_interfaces, n_bins )
    real(dk) :: rates( n_interfaces, n_reactions )
    real(dk) :: expected( n_interfaces, n_reactions )
    real(dk) :: doses( n_interfaces, n_dose_rates )
    real(dk) :: expected_doses( n_interfaces, n_dose_rates )

    config_path = "test/data/solver_from_host.json"
    core => core_t( config_path )

    call assert( 570447568,                                                   &
                 core%number_of_photolysis_reactions( ) == n_reactions )
    call assert( 682765903, core%number_of_dose_rates( ) == n_dose_rates )
    labels = core%photolysis_reaction_labels( )
    call assert( 400290664, labels(1) == "jfoo" )
    call assert( 912659009, labels(2) == "jbar" )

    updater = core%get_radiation_field_updater( found )
    call assert( 742502105, found )

    do i_bin = 1, n_bins
      do i_interface = 1, n_interfaces
        fdr( i_interface, i_bin ) = 0.2_dk * i_interface + 0.05_dk * i_bin
        fup( i_interface, i_bin ) = 0.01_dk * i_interface
        fdn( i_interface, i_bin ) = 0.03_dk * i_bin
        edr( i_interface, i_bin ) = 0.1_dk * i_interface + 0.02_dk * i_bin
        eup( i_interface, i_bin ) = 0.004_dk * i_interface
        edn( i_interface, i_bin ) = 0.006_dk * i_bin
      end do
    end do

    ! TUV-x forms the photolysis rate constants from the actinic flux
    ! components alone, so this call supplies only those three
    call updater%update( direct_actinic_flux   = fdr,                         &
                         upward_actinic_flux   = fup,                         &
                         downward_actinic_flux = fdn )
    call core%run( solar_zenith_angle = 42.0_dk,                              &
                   earth_sun_distance = earth_sun_distance,                   &
                   photolysis_rate_constants = rates,                         &
                   dose_rates = doses )

    ! TUV-x scales the host radiation field by the Earth-Sun distance factor
    ! and multiplies it by the extraterrestrial flux
    call photolysis_reference( fdr, fup, fdn, expected )
    call check_values( 172345201, rates, expected, tol )

    ! TUV-x forms the dose rates from the irradiance components. The host
    ! application omitted them, so TUV-x used zeros and every dose rate is
    ! zero.
    expected_doses(:,:) = 0.0_dk
    call check_values( 553126840, doses, expected_doses, tol )

    ! the reported radiation field carries the Earth-Sun distance factor
    field = core%get_radiation_field( )
    call check_values( 402188297, field%fdr_, fdr * earth_sun_distance, tol )
    call check_values( 914556642, field%fup_, fup * earth_sun_distance, tol )
    call check_values( 744399738, field%fdn_, fdn * earth_sun_distance, tol )

    ! a second call with a new field must give a new answer. This call also
    ! supplies the irradiance components, so the dose rates are non-zero.
    fdr(:,:) = 2.0_dk * fdr(:,:)
    call updater%update( direct_actinic_flux   = fdr,                         &
                         upward_actinic_flux   = fup,                         &
                         downward_actinic_flux = fdn,                         &
                         direct_irradiance     = edr,                         &
                         upward_irradiance     = eup,                         &
                         downward_irradiance   = edn )
    call core%run( solar_zenith_angle = 42.0_dk,                              &
                   earth_sun_distance = earth_sun_distance,                   &
                   photolysis_rate_constants = rates,                         &
                   dose_rates = doses )
    call photolysis_reference( fdr, fup, fdn, expected )
    call check_values( 291557784, rates, expected, tol )

    ! TUV-x converts the normalized irradiance to W m-2 and then applies the
    ! spectral weight of each dose rate
    do i_rate = 1, n_dose_rates
      do i_interface = 1, n_interfaces
        expected_doses( i_interface, i_rate ) = 0.0_dk
        do i_bin = 1, n_bins
          expected_doses( i_interface, i_rate ) =                             &
              expected_doses( i_interface, i_rate ) +                         &
              ( edr( i_interface, i_bin ) + eup( i_interface, i_bin ) +       &
                edn( i_interface, i_bin ) ) * earth_sun_distance * etfl *     &
              hc / ( lambda( i_bin ) * 1.0e-13_dk ) *                         &
              weight( i_bin, i_rate )
        end do
      end do
    end do
    call check_values( 138029475, doses, expected_doses, tol )
    call assert( 585397321, all( doses( :, 1 ) > doses( :, 2 ) ) )

    deallocate( core )

  end subroutine test_core_with_host_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine photolysis_reference( fdr, fup, fdn, expected )
    ! Computes the photolysis rate constants that TUV-x must return for the
    ! given actinic flux components

    real(dk), intent(in)  :: fdr(:,:)      ! direct actinic flux
    real(dk), intent(in)  :: fup(:,:)      ! upward actinic flux
    real(dk), intent(in)  :: fdn(:,:)      ! downward actinic flux
    real(dk), intent(out) :: expected(:,:) ! rate constants

    integer :: i_interface, i_bin, i_reaction

    do i_reaction = 1, n_reactions
      do i_interface = 1, n_interfaces
        expected( i_interface, i_reaction ) = 0.0_dk
        do i_bin = 1, n_bins
          expected( i_interface, i_reaction ) =                               &
              expected( i_interface, i_reaction ) +                           &
              ( fdr( i_interface, i_bin ) + fup( i_interface, i_bin ) +       &
                fdn( i_interface, i_bin ) ) * earth_sun_distance * etfl *     &
              xsqy( i_reaction )
        end do
      end do
    end do

  end subroutine photolysis_reference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_updater_for_another_solver( )
    ! Tests that TUV-x reports a missing updater when the configuration sets
    ! a solver of another type

    use musica_assert,                 only : assert
    use musica_string,                 only : string_t
    use tuvx_core,                     only : core_t
    use tuvx_solver_from_host,         only : radiation_field_updater_t

    class(core_t), pointer :: core
    type(radiation_field_updater_t) :: updater
    type(string_t) :: config_path
    logical :: found

    config_path = "test/data/solver_from_host.other_solver.json"
    core => core_t( config_path )

    updater = core%get_radiation_field_updater( found )
    call assert( 264913507, .not. found )

    deallocate( core )

  end subroutine test_updater_for_another_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine failure_test( code )
    ! Triggers the failure whose assertion code the caller names

    use musica_assert,                 only : die
    use musica_string,                 only : string_t
    use tuvx_core,                     only : core_t
    use tuvx_solver_from_host,         only : radiation_field_updater_t

    character(len=*), intent(in) :: code ! assertion code to trigger

    class(core_t), pointer :: core
    type(radiation_field_updater_t) :: updater
    type(string_t) :: config_path
    logical :: found
    real(dk) :: good( n_interfaces, n_bins )
    real(dk) :: bad( n_interfaces + 1, n_bins )

    good(:,:) = 0.5_dk
    bad(:,:)  = 0.5_dk

    select case( trim( code ) )
    case( "254866372" )
      ! a host-supplied array whose shape does not match the grids
      config_path = "test/data/solver_from_host.json"
      core => core_t( config_path )
      updater = core%get_radiation_field_updater( )
      call updater%update( direct_actinic_flux   = bad,                       &
                           upward_actinic_flux   = good,                      &
                           downward_actinic_flux = good )
    case( "785646610" )
      ! an optional irradiance argument whose shape does not match the grids
      config_path = "test/data/solver_from_host.json"
      core => core_t( config_path )
      updater = core%get_radiation_field_updater( )
      call updater%update( direct_actinic_flux   = good,                      &
                           upward_actinic_flux   = good,                      &
                           downward_actinic_flux = good,                      &
                           direct_irradiance     = bad )
    case( "921177509" )
      ! a request for an updater without the found flag when the
      ! configuration sets a solver of another type
      config_path = "test/data/solver_from_host.other_solver.json"
      core => core_t( config_path )
      updater = core%get_radiation_field_updater( )
    case( "419530284" )
      ! an update through an updater that has no solver behind it
      config_path = "test/data/solver_from_host.other_solver.json"
      core => core_t( config_path )
      updater = core%get_radiation_field_updater( found )
      call updater%update( direct_actinic_flux   = good,                      &
                           upward_actinic_flux   = good,                      &
                           downward_actinic_flux = good )
    case default
      call die( 318628947 )
    end select

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_solver_from_host
