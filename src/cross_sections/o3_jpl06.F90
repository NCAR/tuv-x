! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_o3_jpl06

  use musica_constants,                only : dk => musica_dk
  use tuvx_cross_section,              only : cross_section_t

  implicit none

  private
  public :: cross_section_o3_jpl06_t

  !> Calculator for O3 cross secitons
  !! Calculates O3 cross section using a combination of data sources
  !! and JPL06-recommended temperature dependent values for the
  !! Hartley and Huggins bands
  type, extends(cross_section_t) :: cross_section_o3_jpl06_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_o3_jpl06_t

  interface cross_section_o3_jpl06_t
    module procedure :: constructor
  end interface

  integer, parameter :: FILE_BASE_ID = 1
  integer, parameter :: FILE_LOW_ID  = 2
  integer, parameter :: FILE_HIGH_ID = 3

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Creates an instance of the JPL06 O3 cross section calculator

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: this
    type(config_t),            intent(inout) :: config
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    character(len=*), parameter :: Iam = "JPL06 O3 cross section constructor"
    type(config_t) :: netcdf_file
    type(string_t) :: required_keys(4), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "base netcdf file"
    required_keys(3) = "low-temperature netcdf file"
    required_keys(4) = "high-temperature netcdf file"
    optional_keys(1) = "name"
    call assert_msg( 429850308,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for JPL06 O3 cross "//    &
                     "section." )
    allocate( cross_section_o3_jpl06_t :: this )
    allocate( this%cross_section_parms( 3 ) )

    ! get grid and profile pointers
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ = profile_warehouse%get_ptr( "temperature", "K" )

    ! get base cross section data
    call config%get( "base netcdf file", netcdf_file, Iam )
    call this%process_file( netcdf_file, grid_warehouse,                      &
                            this%cross_section_parms( FILE_BASE_ID ) )

    ! get low-temperature cross section data
    call config%get( "low-temperature netcdf file", netcdf_file, Iam )
    call this%process_file( netcdf_file, grid_warehouse,                      &
                            this%cross_section_parms( FILE_LOW_ID ) )

    ! get high-temperature cross section data
    call config%get( "high-temperature netcdf file", netcdf_file, Iam )
    call this%process_file( netcdf_file, grid_warehouse,                      &
                            this%cross_section_parms( FILE_HIGH_ID ) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculates the cross section for a given set of environmental conditions

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable                     :: cross_section(:,:)
    class(cross_section_o3_jpl06_t), intent(in)    :: this
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse
    logical, optional,               intent(in)    :: at_mid_point

    real(kind=dk), allocatable :: work_cross_section(:,:)
    integer                    :: n_heights, n_wavelengths
    integer                    :: min_wl, max_wl, i_wl, i_height
    class(grid_t), pointer     :: heights, wavelengths
    class(profile_t), pointer  :: temperatures
    real(kind=dk)              :: temperature

    heights => grid_warehouse%get_grid( this%height_grid_ )
    wavelengths => grid_warehouse%get_grid( this%wavelength_grid_ )
    temperatures =>                                                           &
        profile_warehouse%get_profile( this%temperature_profile_ )

    n_heights = heights%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) n_heights = n_heights - 1
    end if
    n_wavelengths = wavelengths%ncells_
    allocate( cross_section( n_heights, n_wavelengths ) )

    min_wl = 1
    do i_wl = 1, wavelengths%ncells_
      if( wavelengths%edge_( i_wl ) .ge. 196.078_dk ) exit
      min_wl = min_wl + 1
    end do
    max_wl = wavelengths%ncells_
    do i_wl = wavelengths%ncells_, 1, -1
      if( wavelengths%edge_( i_wl + 1 ) .le. 342.5_dk ) exit
      max_wl = max_wl - 1
    end do

    do i_height = 1, n_heights
    associate( cs => cross_section( i_height, : ),                            &
               low_cs => this%cross_section_parms( FILE_LOW_ID )%array(:,1),  &
               high_cs => this%cross_section_parms( FILE_HIGH_ID )%array(:,1) )
      temperature = temperatures%edge_val_( i_height )
      if( present( at_mid_point ) ) then
        if( at_mid_point ) temperature = temperatures%mid_val_( i_height )
      end if
      cs(:) = this%cross_section_parms( FILE_BASE_ID )%array(:,1)
      if( temperature .lt. 218.0_dk ) then
        cs( min_wl : max_wl ) = low_cs( min_wl : max_wl )
      else if( temperature .le. 298.0_dk ) then
        cs( min_wl : max_wl ) = low_cs( min_wl : max_wl ) +                   &
            ( high_cs( min_wl : max_wl ) - low_cs( min_wl : max_wl ) )        &
            / ( 298.0_dk - 218.0_dk ) * ( temperature - 218.0_dk )
      else
        cs( min_wl : max_wl ) = high_cs( min_wl : max_wl )
      end if
    end associate
    end do

    deallocate( heights )
    deallocate( wavelengths )
    deallocate( temperatures )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_o3_jpl06
