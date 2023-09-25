! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_no2_tint
! Calculates the cross section for nitrogen dioxide with temperature
! interpolation

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use tuvx_cross_section,              only : cross_section_t

  implicit none

  private
  public :: cross_section_no2_tint_t

  !> Calculator for tint_cross_section
  type, extends(cross_section_t) :: cross_section_no2_tint_t
  contains
    !> Calculate the cross section
    procedure :: calculate
    !> clean up
    final     :: finalize
  end type cross_section_no2_tint_t

  !> Constructor
  interface cross_section_no2_tint_t
    module procedure constructor
  end interface cross_section_no2_tint_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize cross_section_no2_tint_t object

    use musica_assert,                 only : assert_msg, die_msg
    use musica_constants,              only : dk => musica_dk
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(cross_section_no2_tint_t), pointer       :: this ! This :f:type:`~tuvx_cross_section_no2_tint/cross_section_no2_tint_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! local variables
    character(len=*), parameter :: Iam = 'no2 tint cross section constructor'
    character(len=*), parameter :: Hdr = 'cross_section_'

    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk

    integer :: parmNdx, fileNdx, Ndxl, Ndxu
    integer :: nParms, nTemps
    real(dk)              :: tmp
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable :: netcdf_obj
    type(config_t) :: netcdf_files, netcdf_file
    class(iterator_t), pointer :: iter
    type(string_t) :: file_path
    class(grid_t),  pointer     :: lambdaGrid
    type(config_t) :: interpolator_config
    class(interpolator_t), pointer :: interpolator
    type(string_t) :: required_keys(2), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "netcdf files"
    optional_keys(1) = "lower extrapolation"
    optional_keys(2) = "upper extrapolation"
    optional_keys(3) = "name"
    call assert_msg( 309588437,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "no2 tint cross section." )

    allocate(this)

    ! Get model wavelength grids
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ =                                               &
        profile_warehouse%get_ptr( "temperature", "K" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! Get cross section netcdf filespec
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
                                        variable_name = Hdr )
      nParms = size( netcdf_obj%parameters, dim = 2 )
      ! must have at least one parameter
      call assert_msg( 211098593, nParms >= 2,                                &
                     'File: '//file_path//' must have at least 2 parameters.' )
      associate( Xsection => this%cross_section_parms( fileNdx ) )

      call netcdf_file%get( "interpolator", interpolator_config, Iam,         &
                            found = found )
      if( .not. found ) then
        call interpolator_config%empty( )
        call interpolator_config%add( "type", "conserving", Iam )
      end if
      interpolator => interpolator_t( interpolator_config )

      ! interpolation temperatures must be in netcdf file
      call assert_msg( 140564360,  allocated( netcdf_obj%temperature ),       &
                       'File: '//file_path//' does not have interpolation '// &
                       'temperatures.' )
      Xsection%temperature = netcdf_obj%temperature
      nTemps = size( Xsection%temperature )
      ! must have two or more interpolation temperatures
      call assert_msg( 953662959,  nTemps >= 2,                               &
                       'File: '//file_path//' temperature array has less '//  &
                       'than 2 entries.' )
      call assert_msg( 834627068, nTemps >= nParms,                           &
                       'File: '//file_path//' temperature array has less '//  &
                       'than the number of parameters.' )
      Xsection%deltaT = Xsection%temperature( 2 : nParms )                    &
                        - Xsection%temperature( 1 : nParms - 1 )
      monopos = all( Xsection%deltaT > rZERO )
      if( .not. monopos ) then
        call assert_msg( 655847084, .not. any( Xsection%deltaT > rZERO ),     &
                         'File: '//file_path//' temperature array is not '//  &
                         'monotonic.' )
        do Ndxl = 1, nParms / 2
          Ndxu = nParms - Ndxl + 1
          tmp = Xsection%temperature( Ndxl )
          Xsection%temperature( Ndxl ) = Xsection%temperature( Ndxu )
          Xsection%temperature( Ndxu ) = tmp
          data_parameter = netcdf_obj%parameters( :, Ndxl )
          netcdf_obj%parameters( :, Ndxl ) =                                  &
              netcdf_obj%parameters( :, Ndxu )
          netcdf_obj%parameters( :, Ndxu ) = data_parameter
        enddo
        Xsection%deltaT = Xsection%temperature( 2 : nParms )                  &
                          - Xsection%temperature( 1 : nParms - 1 )
      endif

      ! interpolate from data to model wavelength grid
      if( allocated( netcdf_obj%wavelength ) ) then
        if( .not. allocated( Xsection%array ) ) then
          allocate( Xsection%array( lambdaGrid%ncells_, nParms ) )
        endif
        do parmNdx = 1, nParms
          data_lambda    = netcdf_obj%wavelength
          data_parameter = netcdf_obj%parameters( :, parmNdx )
          call this%add_points( config, data_lambda, data_parameter )
          Xsection%array( :, parmNdx ) =                                    &
              interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                        x_source = data_lambda,             &
                                        y_source = data_parameter )
        enddo
      else
        Xsection%array = netcdf_obj%parameters
      endif
      end associate

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

    real(kind=dk), allocatable                     :: cross_section(:,:) ! Calculated cross section
    class(cross_section_no2_tint_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_no2_tint/cross_section_no2_tint_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'no2 tint cross section calculate'
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk

    integer :: nTemp, nzdim
    integer :: fileNdx, tNdx, vertNdx
    real(dk)    :: Tadj, Tstar
    real(dk),         allocatable :: wrkCrossSection(:,:)
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

    allocate( wrkCrossSection( lambdaGrid%ncells_, nzdim ) )
    wrkCrossSection = rZERO

    do fileNdx = 1, size( this%cross_section_parms )
      associate( dataTemp => this%cross_section_parms( fileNdx )%temperature, &
                 wrkXsect => this%cross_section_parms( fileNdx ) )
      nTemp = size( dataTemp )
      do vertNdx = 1, nzdim
        Tadj = min( max( modelTemp( vertNdx ), dataTemp(1) ),                 &
                    dataTemp( nTemp ) )
        do tNdx = 2, nTemp
          if( Tadj <= dataTemp( tNdx ) ) then
            exit
          endif
        enddo
        tNdx = min( nTemp, tNdx ) - 1
        Tstar = ( Tadj - dataTemp( tNdx ) ) / wrkXsect%deltaT( tNdx )
        wrkCrossSection( :, vertNdx ) = wrkCrossSection( :, vertNdx )         &
                                        + wrkXsect%array( :, tNdx )           &
                                    + Tstar * ( wrkXsect%array( :, tNdx + 1 ) &
                                                - wrkXsect%array( :, tNdx ) )
      enddo
      end associate
    enddo

    cross_section = transpose( wrkCrossSection )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! finalize the cross section type

    type(cross_section_no2_tint_t), intent(inout) :: this ! This :f:type:`~tuvx_cross_section_no2_tint/cross_section_no2_tint_t`

    integer :: ndx

    if( allocated( this%cross_section_parms ) ) then
      do ndx = 1, size( this%cross_section_parms )
        associate( Xsection => this%cross_section_parms( ndx ) )
        if( allocated( Xsection%array ) ) then
          deallocate( Xsection%array )
        endif
        if( allocated( Xsection%temperature ) ) then
          deallocate( Xsection%temperature )
        endif
        if( allocated( Xsection%deltaT ) ) then
          deallocate( Xsection%deltaT )
        endif
        end associate
      enddo
      deallocate( this%cross_section_parms )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_no2_tint
