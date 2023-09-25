! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

!> Builds solver_t objects for use in radiative transfer calculations
module tuvx_solver_factory

  use tuvx_solver,                     only : solver_t
  use tuvx_solver_delta_eddington,     only : solver_delta_eddington_t
  use tuvx_solver_discrete_ordinate,   only : solver_discrete_ordinate_t

  implicit none

  private
  public :: solver_builder, solver_type_name, solver_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builder of solver_t objects
  function solver_builder( config, grid_warehouse, profile_warehouse )        &
      result( solver )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    !> Solver configuration
    type(config_t),  intent(inout) :: config
    !> Configured solver
    class(solver_t), pointer       :: solver
    !> Grid warehouse
    type(grid_warehouse_t), intent(in) :: grid_warehouse
    !> Profile warehouse
    type(profile_warehouse_t), intent(in) :: profile_warehouse

    ! Local variables
    type(string_t) :: solver_type
    character(len=*), parameter :: Iam = "solver builder"

    solver => null( )
    call config%get( "type", solver_type, Iam )

    select case( solver_type%to_char( ) )
      case( "delta eddington" )
        solver => solver_delta_eddington_t( config, grid_warehouse,           &
                                            profile_warehouse )
      case( "discrete ordinate" )
        solver => solver_discrete_ordinate_t( config, grid_warehouse,         &
                                              profile_warehouse )
      case default
        call die_msg( 297172205, "Invalid solver type: '"//                   &
                                 solver_type%to_char( ) )
    end select

  end function solver_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the type of a solver as a string
  type(string_t) function solver_type_name( solver ) result( name )

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    !> Solver to return type name for
    class(solver_t), intent(in) :: solver

    select type( solver )
      type is( solver_delta_eddington_t )
        name = "solver_delta_eddington_t"
      type is( solver_discrete_ordinate_t )
        name = "solver_discrete_ordinate_t"
      class default
        call die( 118844516 )
    end select

  end function solver_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a solver pointer as a subclass by type name
  function solver_allocate( type_name ) result( solver )

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    !> Name of the type to allocate
    type(string_t),  intent(in) :: type_name
    !> Allocated, uninitialized  solver
    class(solver_t), pointer    :: solver

    solver => null( )

    select case( type_name%to_char( ) )
      case( "solver_delta_eddington_t" )
        allocate( solver_delta_eddington_t :: solver )
      case( "solver_discrete_ordinate_t" )
        allocate( solver_discrete_ordinate_t :: solver )
      case default
        call die_msg( 947069792, "Invalid solver type name '"//type_name//"'" )
    end select

  end function solver_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_solver_factory
