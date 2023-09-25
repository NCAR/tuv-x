! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiator_factory
! Builds :f:type:`~tuvx_radiator/radiator_t` s for
! :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`.

  use tuvx_radiator,                   only : radiator_t
  use tuvx_radiator_aerosol,           only : radiator_aerosol_t
  use tuvx_radiator_from_host,         only : radiator_from_host_t

  implicit none

  private
  public :: radiator_builder, radiator_type_name, radiator_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function radiator_builder( config, grid_warehouse, profile_warehouse,       &
     cross_section_warehouse ) result( new_radiator )
    ! Builder of :f:type:`~tuvx_radiator/radiator_t` objects

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(config_t),                  intent(inout) :: config        ! Radiator configuration data
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse ! profile warehouse
    type(cross_section_warehouse_t), intent(inout) :: cross_section_warehouse ! cross section warehouse
    class(radiator_t),               pointer       :: new_radiator  ! New :f:type:`~tuvx_radiator/radiator_t` object

    ! Local variables
    type(string_t) :: radiator_type
    character(len=*), parameter :: Iam = 'Radiator builder'

    new_radiator => null()
    call config%get( 'type', radiator_type, Iam )

    select case( radiator_type%to_char() )
      case( 'base' )
        new_radiator => radiator_t( config, grid_warehouse, profile_warehouse,&
                                    cross_section_warehouse )
      case( 'aerosol' )
        new_radiator => radiator_aerosol_t( config, grid_warehouse,           &
                                            profile_warehouse )
      case default
        call die_msg( 460768245, "Invalid radiator type: '"//                 &
                                 radiator_type%to_char()//"'" )
    end select

  end function radiator_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function radiator_type_name( radiator ) result( name )
    ! Returns the type of a radiator as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(radiator_t), intent(in) :: radiator

    select type( radiator )
      type is( radiator_t )
        name = "radiator_t"
      type is( radiator_aerosol_t )
        name = "radiator_aerosol_t"
      type is( radiator_from_host_t )
        name = "radiator_from_host_t"
      class default
        call die( 365718517 )
    end select

  end function radiator_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function radiator_allocate( type_name ) result( radiator )
    ! Allocates a radiator pointer as the base or a subclass by type name

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t), intent(in) :: type_name ! Name of the type to allocate
    class(radiator_t), pointer :: radiator  ! Allocated radiator

    radiator => null( )

    select case( type_name%to_char( ) )
      case( 'radiator_t' )
        allocate( radiator_t :: radiator )
      case( 'radiator_aerosol_t' )
        allocate( radiator_aerosol_t :: radiator )
      case( 'radiator_from_host_t' )
        allocate( radiator_from_host_t :: radiator )
      case default
        call die_msg( 670539061, "Invalid radiator type: '"//                 &
                                 type_name//"'" )
    end select

  end function radiator_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_factory
