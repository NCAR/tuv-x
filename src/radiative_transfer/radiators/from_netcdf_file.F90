! Copyright (C) 2023 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module tuvx_radiator_from_netcdf_file

  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t
  use tuvx_radiator,                   only : radiator_t

  implicit none

  private
  public :: radiator_from_netcdf_file_t

  !> User-specified radiator
  !!
  !! User-specified radiators have fixed optical properties that are read
  !! in from a NetCDF file.
  type, extends(radiator_t) :: radiator_from_netcdf_file_t
  contains
    !> Update the radiator for new environmental conditions
    procedure :: update_state
  end type radiator_from_netcdf_file_t

  !> User-specified radiator constructor
  interface radiator_from_netcdf_file_t
    module procedure constructor
  end interface radiator_from_netcdf_file_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a user-specified radiator
  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use musica_io_netcdf,              only : io_netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(radiator_from_netcdf_file_t), pointer       :: this ! The constructor radiator
    type(config_t),                     intent(inout) :: config ! configuration for the radiator
    type(grid_warehouse_t),             intent(inout) :: grid_warehouse ! Configured grids
    type(profile_warehouse_t),          intent(inout) :: profile_warehouse ! Configured profiles

    character(len=*),  parameter :: Iam =                                     &
                                    "radiator from NetCDF file constructor"
    class(grid_t), pointer :: heights, wavelengths
    type(string_t) :: file_path, variable
    type(io_netcdf_t), pointer :: netcdf_file
    real(kind=dk), allocatable :: temp_G(:,:)
    integer :: i_wl, i_height
    type(string_t) :: required_keys(3), optional_keys(0)

    required_keys(1) = "type"
    required_keys(2) = "name"
    required_keys(3) = "file path"
    call assert_msg( 723245326,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for radiator from NetCDF" )

    allocate(this)

    call config%get( "name", this%handle_, Iam )
    call config%get( "type", this%type_,   Iam )

    heights => grid_warehouse%get_grid( "height", "km" )
    wavelengths => grid_warehouse%get_grid( "wavelength", "nm" )

    allocate( this%state_%layer_OD_(  heights%ncells_, wavelengths%ncells_ ) )
    allocate( this%state_%layer_SSA_( heights%ncells_, wavelengths%ncells_ ) )
    allocate( this%state_%layer_G_( heights%ncells_, wavelengths%ncells_, 1 ) )
    allocate( temp_G( heights%ncells_, wavelengths%ncells_ ) )

    call config%get( "file path", file_path, Iam )
    netcdf_file => io_netcdf_t( file_path, read_only = .true. )
    variable = "optical_depth"
    call netcdf_file%read_2D_double( variable, this%state_%layer_OD_,  Iam )
    variable = "single_scattering_albedo"
    call netcdf_file%read_2D_double( variable, this%state_%layer_SSA_, Iam )
    variable = "asymmetry_factor"
    call netcdf_file%read_2D_double( variable, temp_G, Iam )
    this%state_%layer_G_(:,:,1) = temp_G(:,:)

    deallocate( netcdf_file )
    deallocate( heights )
    deallocate( wavelengths )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the radiator state for current conditions
  subroutine update_state( this, grid_warehouse, profile_warehouse,           &
      cross_section_warehouse )

    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(radiator_from_netcdf_file_t), intent(inout) :: this ! Radiator to update
    type(grid_warehouse_t),             intent(inout) :: grid_warehouse ! Current grids
    type(profile_warehouse_t),          intent(inout) :: profile_warehouse ! Current profiles
    type(cross_section_warehouse_t),    intent(inout) :: cross_section_warehouse ! Current cross-sections

    ! nothing to update - radiator properties for this type are fixed

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_from_netcdf_file