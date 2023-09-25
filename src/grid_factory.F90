! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_factory
! Provides a function which creates grids for the
! :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`.

  use tuvx_grid,                       only : grid_t
  use tuvx_grid_equal_delta,           only : grid_equal_delta_t
  use tuvx_grid_from_csv_file,         only : grid_from_csv_file_t
  use tuvx_grid_from_config,           only : grid_from_config_t
  use tuvx_grid_from_host,             only : grid_from_host_t

  implicit none

  private
  public :: grid_builder, grid_type_name, grid_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function grid_builder( config ) result( new_grid_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> Arguments
    !> Grid configuration data
    type(config_t), intent(inout) :: config

    class(grid_t), pointer :: new_grid_t ! New A :f:type:`~tuvx_grid/grid_t` object

    !> Local variables
    type(string_t) :: grid_type
    character(len=*), parameter :: Iam = 'Grid builder: '

    new_grid_t => null()
    call config%get( 'type', grid_type, Iam )

    select case( grid_type%to_char() )
      case( 'equal interval' )
        new_grid_t => grid_equal_delta_t( config )
      case( 'from csv file' )
        new_grid_t => grid_from_csv_file_t( config )
      case( 'from config file' )
        new_grid_t => grid_from_config_t( config )
      case default
        call die_msg( 460768215, "Invalid grid type: '" &
          // grid_type%to_char()//"'" )
    end select

  end function grid_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function grid_type_name( grid ) result( name )
    ! Returns the type of a grid as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(grid_t), intent(in) :: grid

    select type( grid )
      type is( grid_equal_delta_t )
        name = "grid_equal_delta_t"
      type is( grid_from_csv_file_t )
        name = "grid_from_csv_file_t"
      type is( grid_from_config_t )
        name = "grid_from_config_t"
      type is( grid_from_host_t )
        name = "grid_from_host_t"
      class default
        call die( 983843127 )
    end select

  end function grid_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function grid_allocate( type_name ) result( grid )
    ! Allocates a grid pointer as a subclass by type

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t), intent(in) :: type_name ! name of the type to allocate
    class(grid_t),  pointer    :: grid      ! allocated grid

    grid => null( )

    select case( type_name%to_char( ) )
      case( 'grid_equal_delta_t' )
        allocate( grid_equal_delta_t :: grid )
      case( 'grid_from_csv_file_t' )
        allocate( grid_from_csv_file_t :: grid )
      case( 'grid_from_config_t' )
        allocate( grid_from_config_t :: grid )
      case( 'grid_from_host_t' )
        allocate( grid_from_host_t :: grid )
      case default
        call die_msg( 351430046, "Invalid grid type: '"//type_name//"'" )
    end select

  end function grid_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_factory
