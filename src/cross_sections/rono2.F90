! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_rono2
! Calculates the cross section for a nitrate ester

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_cross_section,              only : cross_section_t

  implicit none

  private
  public :: cross_section_rono2_t

  !> Calculator for rono2 cross section
  type, extends(cross_section_t) :: cross_section_rono2_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_rono2_t

  !> Constructor
  interface cross_section_rono2_t
    module procedure constructor
  end interface cross_section_rono2_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize cross_section_t object

    use musica_assert,                 only : assert_msg, die_msg
    use musica_constants,              only : dk => musica_dk
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_rono2_t),    pointer :: this ! This :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'rono2 cross section initialize'
    character(len=*), parameter :: Hdr = 'cross_section_'
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk

    integer :: parmNdx, fileNdx, nParms
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found
    character(len=:), allocatable :: msg
    character(len=:), allocatable :: addpntKey
    type(netcdf_t),   allocatable :: netcdf_obj
    type(config_t) :: netcdf_files, netcdf_file
    class(iterator_t), pointer :: iter
    type(string_t) :: file_path
    type(config_t)                :: tmp_config, extrap_config
    class(grid_t),    pointer     :: lambdaGrid
    type(config_t) :: interpolator_config
    class(interpolator_t), pointer :: interpolator
    type(string_t) :: required_keys(2), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "netcdf files"
    optional_keys(1) = "lower extrapolation"
    optional_keys(2) = "upper extrapolation"
    optional_keys(3) = "name"
    call assert_msg( 191004841,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "rono2 cross section." )

    allocate( this )

    ! Get model wavelength grids
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ =                                               &
        profile_warehouse%get_ptr( "temperature", "K" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! get cross section netcdf filespec
    call config%get( 'netcdf files', netcdf_files, Iam )
    iter => netcdf_files%get_iterator( )

    allocate( this%cross_section_parms( netcdf_files%number_of_children( ) ) )

    fileNdx = 0
    do while( iter%next( ) )
      call netcdf_files%get( iter, netcdf_file, Iam )
      fileNdx = fileNdx + 1

      allocate( netcdf_obj )
      ! read netcdf cross section parameters
      call netcdf_file%get( "file path", file_path, Iam )
      call netcdf_obj%read_netcdf_file( file_path = file_path%to_char( ),     &
                                        variable_name  = Hdr )
      nParms = size( netcdf_obj%parameters, dim = 2 )
      call assert_msg( 236151624, nParms >= 1,                                &
                       'File: '//file_path//' contains no parameters.' )

      call netcdf_file%get( "interpolator", interpolator_config, Iam,         &
                            found = found )
      if( .not. found ) then
        call interpolator_config%empty( )
        call interpolator_config%add( "type", "conserving", Iam )
      end if
      interpolator => interpolator_t( interpolator_config )

      ! interpolate from data to model wavelength grid
      if( allocated( netcdf_obj%wavelength ) ) then
        if( .not. allocated( this%cross_section_parms( fileNdx )%array ) )    &
            then
          allocate( this%cross_section_parms(                                 &
                             fileNdx )%array( lambdaGrid%ncells_, nParms ) )
        endif
        do parmNdx = 1, nParms
          data_lambda    = netcdf_obj%wavelength
          data_parameter = netcdf_obj%parameters( :, parmNdx )
          if( parmNdx == 1 ) then
            call this%add_points( config, data_lambda, data_parameter )
          elseif( parmNdx == 2 ) then
            tmp_config = config
            addpntKey = 'lower extrapolation'
            call extrap_config%empty( )
            call extrap_config%add( 'type', 'boundary', Iam )
            call tmp_config%add( addpntKey, extrap_config, Iam )
            addpntKey = 'upper extrapolation'
            call tmp_config%add( addpntKey, extrap_config, Iam )
            call this%add_points( tmp_config, data_lambda, data_parameter )
          endif
          this%cross_section_parms( fileNdx )%array( :, parmNdx ) =           &
              interpolator%interpolate( x_target = lambdaGrid%edge_,          &
                                        x_source = data_lambda,               &
                                        y_source = data_parameter )
        enddo
      else
        this%cross_section_parms( fileNdx )%array = netcdf_obj%parameters
      endif
      if( allocated( netcdf_obj%temperature ) ) then
        this%cross_section_parms( fileNdx )%temperature =                     &
            netcdf_obj%temperature
      endif
      deallocate( interpolator )
      deallocate( netcdf_obj   )
    enddo

    deallocate( iter       )
    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable                  :: cross_section(:,:) ! Calculated cross section
    class(cross_section_rono2_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_rono2/cross_section_rono2_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam = 'rono2 cross section calculate'
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk

    integer :: vertNdx
    integer :: nzdim
    real(dk), parameter :: T0 = 298._dk
    real(dk) :: Temp
    real(dk),         allocatable :: modelTemp(:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
        modelTemp = mdlTemperature%mid_val_
      else
        modelTemp = mdlTemperature%edge_val_
      endif
    else
      modelTemp = mdlTemperature%edge_val_
    endif

    allocate( cross_section( lambdaGrid%ncells_, nzdim ) )
    cross_section = rZERO

    do vertNdx = 1, nzdim
      Temp = modelTemp( vertNdx ) - T0
      cross_section( :, vertNdx ) = this%cross_section_parms(1)%array(:,1)    &
                        * exp( this%cross_section_parms(1)%array(:,2) * Temp )
    enddo

    cross_section = transpose( cross_section )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_rono2
